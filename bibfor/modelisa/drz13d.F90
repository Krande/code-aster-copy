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

subroutine drz13d(noma, ligrmo, type_vale, nb_node, list_node, &
                  cmp_index_dx, cmp_index_dy, cmp_index_dz, cmp_index_drx, cmp_index_dry, &
                  cmp_index_drz, lisrel, nom_noeuds)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/afrela.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisdg.h"
#include "asterfort/int_to_char8.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
!
!
    character(len=8), intent(in) :: noma
    character(len=19), intent(in) :: ligrmo
    character(len=4), intent(in) :: type_vale
    integer(kind=8), intent(in) :: nb_node
    character(len=24), intent(in) :: list_node
    integer(kind=8), intent(in) :: cmp_index_dx
    integer(kind=8), intent(in) :: cmp_index_dy
    integer(kind=8), intent(in) :: cmp_index_dz
    integer(kind=8), intent(in) :: cmp_index_drx
    integer(kind=8), intent(in) :: cmp_index_dry
    integer(kind=8), intent(in) :: cmp_index_drz
    character(len=19), intent(in) :: lisrel
    character(len=8), intent(out) :: nom_noeuds(:)
!
! --------------------------------------------------------------------------------------------------
!
! Loads - Affectation
!
! Apply transformation - 3D with at least one node has DRX, DRY and DRZ dof
!
! --------------------------------------------------------------------------------------------------
!
! In  noma          : mesh
! In  ligrmo        : <LIGREL> of model
! In  type_vale     : type of affected value
! In  nb_node       : number of nodes  applying translation
! In  list_node     : list of nodes applying translation
! In  tran          : vector defining translation
! In  cmp_index_dx  : index in DEPL_R <GRANDEUR> for DX
! In  cmp_index_dy  : index in DEPL_R <GRANDEUR> for DY
! In  cmp_index_dz  : index in DEPL_R <GRANDEUR> for DZ
! In  cmp_index_drx : index in DEPL_R <GRANDEUR> for DRX
! In  cmp_index_dry : index in DEPL_R <GRANDEUR> for DRY
! In  cmp_index_drz : index in DEPL_R <GRANDEUR> for DRZ
! In  lisrel        : list of relations
! Out nom_noeuds    : nom des noeuds "maitres" pour la relation
!
! --------------------------------------------------------------------------------------------------
!
!
    integer(kind=8) :: i_no
    integer(kind=8) ::   jprnm
    integer(kind=8) :: nbec
    integer(kind=8) :: jlino, numnoe_m, numnoe_a
    integer(kind=8) :: nb_maxi, nb_term, ier
    real(kind=8) :: un, x, y, z
    real(kind=8) :: vale_real
    complex(kind=8) :: vale_cplx
    character(len=8) :: vale_fonc
    character(len=4) :: type_coef
    character(len=8) :: nomg, nomnoe_m, nomnoe_a
    complex(kind=8), pointer :: coec(:) => null()
    real(kind=8), pointer :: coer(:) => null()
    integer(kind=8), pointer :: dime(:) => null()
    real(kind=8), pointer :: direct(:) => null()
    character(len=8), pointer :: lisddl(:) => null()
    character(len=8), pointer :: lisno(:) => null()
    real(kind=8), pointer :: vale(:) => null()
    aster_logical :: lcolle
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Initializations
!
    vale_fonc = '&FOZERO'
    vale_real = 0.d0
    vale_cplx = (0.d0, 0.d0)
    un = 1.d0
    type_coef = 'REEL'
    ASSERT(type_vale .ne. 'COMP')
!
! - Information about <GRANDEUR>
!
    nomg = 'DEPL_R'
    call dismoi('NB_EC', nomg, 'GRANDEUR', repi=nbec)
    ASSERT(nbec .le. 11)
    call jeveuo(ligrmo//'.PRNM', 'L', jprnm)
!
! - Nodes coordinates
!
    call jeveuo(noma//'.COORDO    .VALE', 'L', vr=vale)
!
! - List of nodes to apply linear relation
!
    call jeveuo(list_node, 'L', jlino)
!
! - Working vectors
!
    nb_maxi = 4
    AS_ALLOCATE(vk8=lisno, size=nb_maxi)
    AS_ALLOCATE(vk8=lisddl, size=nb_maxi)
    AS_ALLOCATE(vr=coer, size=nb_maxi)
    AS_ALLOCATE(vc=coec, size=nb_maxi)
    AS_ALLOCATE(vr=direct, size=3*nb_maxi)
    AS_ALLOCATE(vi=dime, size=nb_maxi)
!
! - First node with DRX/DRY/DRZ node (reference)
!
    do i_no = 1, nb_node
        numnoe_m = zi(jlino+i_no-1)
        if (exisdg(zi(jprnm-1+(numnoe_m-1)*nbec+1), cmp_index_drx) .and. &
            exisdg(zi(jprnm-1+(numnoe_m-1)*nbec+1), cmp_index_dry) .and. &
            exisdg(zi(jprnm-1+(numnoe_m-1)*nbec+1), cmp_index_drz)) then
            numnoe_a = numnoe_m
            goto 30
        end if
    end do
!
! - No node with DRX/DRY/DRZ: IMPOSSIBLE !
!
    ASSERT(.false.)
!
30  continue
!
    lcolle = .false.
    call jeexin(noma//'.NOMNOE', ier)
    if (ier .ne. 0) then
        lcolle = .true.
    end if
    nomnoe_a = int_to_char8(numnoe_a, lcolle, noma, 'NOEUD')
    nom_noeuds(1) = nomnoe_a
!
! - Loop on nodes
!
    do i_no = 1, nb_node
        numnoe_m = zi(jlino+i_no-1)
        nomnoe_m = int_to_char8(numnoe_m, lcolle, noma, 'NOEUD')
!
        if (numnoe_m .ne. numnoe_a) then
!
! --------- Distances: x = DX(A) - DX(M) and y = DY(A) - DY(M) and z = DZ(A) - DZ(M)
!
            x = vale(3*(numnoe_m-1)+1)-vale(3*(numnoe_a-1)+1)
            y = vale(3*(numnoe_m-1)+2)-vale(3*(numnoe_a-1)+2)
            z = vale(3*(numnoe_m-1)+3)-vale(3*(numnoe_a-1)+3)
!
! --------- Linear relations for translation dof
!
            if (exisdg(zi(jprnm-1+(numnoe_m-1)*nbec+1), cmp_index_dx) .and. &
                exisdg(zi(jprnm-1+(numnoe_m-1)*nbec+1), cmp_index_dy) .and. &
                exisdg(zi(jprnm-1+(numnoe_m-1)*nbec+1), cmp_index_dz)) then
!
                nb_term = 4
                lisno(1) = nomnoe_m
                lisno(2) = nomnoe_a
                lisno(3) = nomnoe_a
                lisno(4) = nomnoe_a
!
! ------------- First relation: DX(M) - DX(A) - Z*DRY(A) + Y*DRZ(A) =0
!
                lisddl(1) = 'DX'
                lisddl(2) = 'DX'
                lisddl(3) = 'DRY'
                lisddl(4) = 'DRZ'
                coer(1) = un
                coer(2) = -un
                coer(3) = -z
                coer(4) = y
!
! ------------- Compute linear relation
!
                call afrela(coer, coec, lisddl, lisno, dime, &
                            direct, nb_term, vale_real, vale_cplx, vale_fonc, &
                            type_coef, type_vale, 0.d0, lisrel)
!
! ------------- Second relation: DY(M) - DY(A) - X*DRZ(A) + Z*DRX(A) =0
!
                lisddl(1) = 'DY'
                lisddl(2) = 'DY'
                lisddl(3) = 'DRZ'
                lisddl(4) = 'DRX'
                coer(1) = un
                coer(2) = -un
                coer(3) = -x
                coer(4) = z
!
! ------------- Compute linear relation
!
                call afrela(coer, coec, lisddl, lisno, dime, &
                            direct, nb_term, vale_real, vale_cplx, vale_fonc, &
                            type_coef, type_vale, 0.d0, lisrel)
!
! ------------- Third relation: DZ(M) - DZ(A) - Y*DRX(A) + X*DRY(A) =0
!
                lisddl(1) = 'DZ'
                lisddl(2) = 'DZ'
                lisddl(3) = 'DRX'
                lisddl(4) = 'DRY'
                coer(1) = un
                coer(2) = -un
                coer(3) = -y
                coer(4) = x
!
! ------------- Compute linear relation
!
                call afrela(coer, coec, lisddl, lisno, dime, &
                            direct, nb_term, vale_real, vale_cplx, vale_fonc, &
                            type_coef, type_vale, 0.d0, lisrel)
            end if
!
! --------- Linear relations for rotation dof
!
            if (exisdg(zi(jprnm-1+(numnoe_m-1)*nbec+1), cmp_index_drx) .and. &
                exisdg(zi(jprnm-1+(numnoe_m-1)*nbec+1), cmp_index_dry) .and. &
                exisdg(zi(jprnm-1+(numnoe_m-1)*nbec+1), cmp_index_drz)) then
!
                nb_term = 2
                lisno(1) = nomnoe_m
                lisno(2) = nomnoe_a
!
! ------------- Fourth relation: DRX(M) - DRX(A)  = 0
!
                lisddl(1) = 'DRX'
                lisddl(2) = 'DRX'
                coer(1) = un
                coer(2) = -un
!
! ------------- Compute linear relation
!
                call afrela(coer, coec, lisddl, lisno, dime, &
                            direct, nb_term, vale_real, vale_cplx, vale_fonc, &
                            type_coef, type_vale, 0.d0, lisrel)
!
!
! ------------- Fifth relation: DRY(M) - DRY(A)  = 0
!
                lisddl(1) = 'DRY'
                lisddl(2) = 'DRY'
                coer(1) = un
                coer(2) = -un
!
! ------------- Compute linear relation
!
                call afrela(coer, coec, lisddl, lisno, dime, &
                            direct, nb_term, vale_real, vale_cplx, vale_fonc, &
                            type_coef, type_vale, 0.d0, lisrel)
!
! ------------- Sixth relation: DRZ(M) - DRZ(A)  = 0
!
                lisddl(1) = 'DRZ'
                lisddl(2) = 'DRZ'
                coer(1) = un
                coer(2) = -un
!
! ------------- Compute linear relation
!
                call afrela(coer, coec, lisddl, lisno, dime, &
                            direct, nb_term, vale_real, vale_cplx, vale_fonc, &
                            type_coef, type_vale, 0.d0, lisrel)
            end if
        end if
    end do
!
    AS_DEALLOCATE(vk8=lisno)
    AS_DEALLOCATE(vk8=lisddl)
    AS_DEALLOCATE(vr=coer)
    AS_DEALLOCATE(vc=coec)
    AS_DEALLOCATE(vr=direct)
    AS_DEALLOCATE(vi=dime)
!
    call jedema()
end subroutine
