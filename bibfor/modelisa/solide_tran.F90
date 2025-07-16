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

subroutine solide_tran(type_geo, noma, type_vale, dist_mini, nb_node, list_node, &
                       lisrel, nom_noeuds, dim)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/afrela.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/getcara_lisno.h"
#include "asterfort/coor_bary.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/int_to_char8.h"
!
!
    character(len=2), intent(in)  :: type_geo
    character(len=8), intent(in)  :: noma
    character(len=4), intent(in)  :: type_vale
    real(kind=8), intent(in)      :: dist_mini
    integer(kind=8), intent(in)           :: nb_node
    character(len=24), intent(in) :: list_node
    character(len=19), intent(in) :: lisrel
    character(len=8), intent(out) :: nom_noeuds(:)
    integer(kind=8), intent(out)          :: dim
!
! --------------------------------------------------------------------------------------------------
!
! Loads - Affectation
!
! Apply transformation - 3D with without nodes with DRZ dof
!
! --------------------------------------------------------------------------------------------------
!
! In  type_geo      : '2D' / '3D'
! In  noma          : mesh
! In  type_vale     : type of affected value
! In  dist_mini     : minimum distance to detect nodes in same place
! In  nb_node       : number of nodes  applying translation
! In  list_node     : list of nodes applying translation
! In  lisrel        : list of relations
! Out nom_noeuds    : nom des (dim+1) noeuds "maitres"
! Out dim           : "dimension" du solide : 0/1/2/3
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: k, km, ka, kb
    integer(kind=8) :: numnoe_m, numnoe_a, numnoe_b
    character(len=8) :: nomnoe_m, nomnoe_a, nomnoe_b
    integer(kind=8) ::    jlino, ier

    integer(kind=8) :: nb_maxi, nb_term, linocara(4), nbnot
    real(kind=8) :: un, cobary(4)
    real(kind=8) :: xa, ya, xb, yb, za, zb
    real(kind=8) :: vale_real
    complex(kind=8) :: vale_cplx
    character(len=8) :: vale_fonc
    real(kind=8) :: xm(3)
    character(len=4) :: type_coef
    aster_logical :: l3d, lcolle

    complex(kind=8), pointer :: coec(:) => null()
    real(kind=8), pointer :: coer(:) => null()
    integer(kind=8), pointer :: dime(:) => null()
    real(kind=8), pointer :: direct(:) => null()
    character(len=8), pointer :: lisddl(:) => null()
    character(len=8), pointer :: lisno(:) => null()
    real(kind=8), pointer :: coor(:) => null()
    integer(kind=8), pointer :: nunocara(:) => null()
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
    ASSERT(dist_mini .gt. 0.d0)
    ASSERT(type_geo .eq. '2D' .or. type_geo .eq. '3D')
    l3d = type_geo .eq. '3D'
!
! - Nodes coordinates
!
    call jeveuo(noma//'.COORDO    .VALE', 'L', vr=coor)
    nbnot = size(coor)/3
    lcolle = .false.
    call jeexin(noma//'.NOMNOE', ier)
    if (ier .ne. 0) then
        lcolle = .true.
    end if
!
! - List of nodes to apply linear relation
!
    call jeveuo(list_node, 'L', jlino)
!
! - Working vectors
!
    nb_maxi = 10
    AS_ALLOCATE(vk8=lisno, size=nb_maxi)
    AS_ALLOCATE(vk8=lisddl, size=nb_maxi)
    AS_ALLOCATE(vr=coer, size=nb_maxi)
    AS_ALLOCATE(vc=coec, size=nb_maxi)
    AS_ALLOCATE(vr=direct, size=3*nb_maxi)
    AS_ALLOCATE(vi=dime, size=nb_maxi)

!   -- Quelle est la situation geometrique ?
!   -----------------------------------------
    call getcara_lisno(noma, nb_node, zi(jlino), dist_mini, dim, linocara)
    ASSERT(dim .le. 3)

!   -- 1) Les relations potentiellement non-lineaires sont celles traduisant
!         l'indeformabilite des dim+1 noeuds de linocara
!   ----------------------------------------------------------------------

!   -- boucle sur les couples de points A, B de linocara :
!   ---------------------------------------------------------
    do ka = 1, dim
        numnoe_a = linocara(ka)
        nomnoe_a = int_to_char8(numnoe_a, lcolle, noma, 'NOEUD')
        xa = coor(3*(numnoe_a-1)+1)
        ya = coor(3*(numnoe_a-1)+2)
        if (l3d) za = coor(3*(numnoe_a-1)+3)

        do kb = ka+1, dim+1
            numnoe_b = linocara(kb)
            nomnoe_b = int_to_char8(numnoe_b, lcolle, noma, 'NOEUD')
            xb = coor(3*(numnoe_b-1)+1)
            yb = coor(3*(numnoe_b-1)+2)
            if (l3d) zb = coor(3*(numnoe_b-1)+3)

            if (l3d) then
                nb_term = 6
            else
                nb_term = 4
            end if

!           -- Relation: AB^2 = cste

!           -- Ordre : A,   B,     A,   B      A,   B
!                     'DX','DX',  'DY','DY',  'DZ','DZ'

            lisno(1) = nomnoe_a
            lisno(2) = nomnoe_b
            lisno(3) = nomnoe_a
            lisno(4) = nomnoe_b
            if (l3d) then
                lisno(5) = nomnoe_a
                lisno(6) = nomnoe_b
            end if

            lisddl(1) = 'DX'
            lisddl(2) = 'DX'
            lisddl(3) = 'DY'
            lisddl(4) = 'DY'
            if (l3d) then
                lisddl(5) = 'DZ'
                lisddl(6) = 'DZ'
            end if

            coer(1) = -2*(xb-xa)
            coer(2) = 2*(xb-xa)
            coer(3) = -2*(yb-ya)
            coer(4) = 2*(yb-ya)
            if (l3d) then
                coer(5) = -2*(zb-za)
                coer(6) = 2*(zb-za)
            end if
!
! --------- Add new linear relation
! --------- Warning epsi=-1.d0 to keep ALL coefficients even there are zero ! (see issue23299)
!
            call afrela(coer, coec, lisddl, lisno, dime, &
                        direct, nb_term, vale_real, vale_cplx, vale_fonc, &
                        type_coef, type_vale, -1.d0, lisrel)
        end do

    end do
    if (nb_node .eq. dim+1) goto 999

!   -- 2) Les relations restantes sont toujours lineaires.
!   ----------------------------------------------------------------------

!   -- boucle sur les noeuds M n'appartenant pas a linocara :
!      On exprime que M est relie aux dim+1 noeuds de linocara
!   ---------------------------------------------------------
    nb_term = dim+2

    AS_ALLOCATE(vi=nunocara, size=nbnot)
    nunocara = 0
    do k = 1, dim+1
        numnoe_a = linocara(k)
        nomnoe_a = int_to_char8(numnoe_a, lcolle, noma, 'NOEUD')
        lisno(1+k) = nomnoe_a
        nunocara(numnoe_a) = 1
    end do

    do km = 1, nb_node
        numnoe_m = zi(jlino-1+km)
        if (nunocara(numnoe_m) .eq. 1) cycle

        xm(1:3) = coor(3*(numnoe_m-1)+1:3*(numnoe_m-1)+3)

        call coor_bary(coor, xm, dim, linocara, cobary)

        nomnoe_m = int_to_char8(numnoe_m, lcolle, noma, 'NOEUD')
        lisno(1) = nomnoe_m

        coer(1) = -1.d0
        do k = 1, dim+1
            coer(1+k) = cobary(k)
        end do

!       -- relation pour DX :
        lisddl(1:dim+2) = 'DX'
        call afrela(coer, coec, lisddl, lisno, dime, &
                    direct, nb_term, vale_real, vale_cplx, vale_fonc, &
                    type_coef, type_vale, 0.d0, lisrel)

!       -- relation pour DY :
        lisddl(1:dim+2) = 'DY'
        call afrela(coer, coec, lisddl, lisno, dime, &
                    direct, nb_term, vale_real, vale_cplx, vale_fonc, &
                    type_coef, type_vale, 0.d0, lisrel)

!       -- relation pour DZ :
        if (l3d) then
            lisddl(1:dim+2) = 'DZ'
            call afrela(coer, coec, lisddl, lisno, dime, &
                        direct, nb_term, vale_real, vale_cplx, vale_fonc, &
                        type_coef, type_vale, 0.d0, lisrel)
        end if

    end do

999 continue

!   -- remplissage de nom_noeuds :
!   ----------------------------------------------
    ASSERT(size(nom_noeuds) .ge. dim+1)
    do k = 1, dim+1
        nomnoe_a = int_to_char8(linocara(k), lcolle, noma, 'NOEUD')
        nom_noeuds(k) = nomnoe_a
    end do

    AS_DEALLOCATE(vi=nunocara)
    AS_DEALLOCATE(vk8=lisno)
    AS_DEALLOCATE(vk8=lisddl)
    AS_DEALLOCATE(vr=coer)
    AS_DEALLOCATE(vc=coec)
    AS_DEALLOCATE(vr=direct)
    AS_DEALLOCATE(vi=dime)

    call jedema()
end subroutine
