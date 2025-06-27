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

subroutine aceagb(nomu, noma, locamb, nbocc)

    use tenseur_dime_module, only: prod_vect
    implicit none
#include "asterf_types.h"
#include "asterc/r8rddg.h"
#include "asterc/r8dgrd.h"
#include "asterfort/alcart.h"
#include "asterfort/angvx.h"
#include "asterfort/assert.h"
#include "asterfort/exisd.h"
#include "asterfort/getvem.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/nocart.h"
#include "asterfort/utmess.h"

!
    integer(kind=8), intent(in) :: nbocc
    aster_logical, intent(in) :: locamb
    character(len=8), intent(in) :: nomu, noma
! --------------------------------------------------------------------------------------------------
!                          AFFE_CARA_ELEM
!
!     AFFECTATION DES CARACTERISTIQUES POUR LE MOT CLE "GRILLE"
!
! --------------------------------------------------------------------------------------------------
!  IN
!     nomu   : nom utilisateur de la commande
!     noma   : nom du maillage
!     locamb : si éléments membrane dans la modelisation
!     nbocc  : nombre d'occurences du mot cle GRILLE
! --------------------------------------------------------------------------------------------------
    real(kind=8), parameter:: epsi = 1.d-6
! --------------------------------------------------------------------------------------------------
    integer(kind=8) :: ioc, nbr_gma, ibid, n, iret, ima, numa, igr, nno, noe
    integer(kind=8) :: mc_section, mc_angl_rep_1, mc_angl_rep_2
    integer(kind=8) :: mc_excentrement, mc_vect_1, mc_vect_2
    integer(kind=8) :: mc_orig, mc_axe_z, mc_section_fo, mc_excentrement_fo
    real(kind=8) :: ang1(2), ang2(2), angrd(2), sl, ez, axex(3), axey(3), norm
    real(kind=8) :: vect_er(3), vect_et(3), vect_ez(3), orig(3)
    real(kind=8) :: vtg1(3), vtg2(3), vec_n(3), bary(3), coor(3, 4), proj(3), vec_x(3), vec_y(3)
    character(len=8) :: slf, ezf
    character(len=19) :: cartgr, cartcf
    character(len=24) :: repere, k24bid
    aster_logical :: lcartf, tria, quad
    character(len=24), allocatable:: l_gma(:)
    integer(kind=8), pointer :: gma(:) => null()
    integer(kind=8), pointer :: l_nds(:) => null()
    real(kind=8), pointer :: val_cmp_r(:) => null()
    real(kind=8), pointer :: coor_all(:) => null()
    character(len=8), pointer :: nom_cmp_r(:) => null()
    character(len=8), pointer :: nom_cmp_f(:) => null()
    character(len=8), pointer :: val_cmp_f(:) => null()
! --------------------------------------------------------------------------------------------------

    call jemarq()

    repere = 'GLOBAL'

    ! Coordonnées des noeuds
    call jeveuo(noma//'.COORDO    .VALE', 'L', vr=coor_all)

! --------------------------------------------------------------------------------------------------
!  ALLOCATION DES CARTES
! --------------------------------------------------------------------------------------------------

    ! Carte pour les valeurs réelles
    cartgr = nomu//'.CARCOQUE'
    call exisd('CARTE', cartgr, iret)
    if (iret .eq. 0) call alcart('G', cartgr, noma, 'CACOQU_R')

    ! Objets intermédiaires pour l'allocation de valeurs
    call jeveuo(cartgr//'.NCMP', 'E', vk8=nom_cmp_r)
    call jeveuo(cartgr//'.VALV', 'E', vr=val_cmp_r)
    ASSERT(size(nom_cmp_r) .ge. 5)
    nom_cmp_r(1:5) = ['SECT_L  ', 'ALPHA   ', 'BETA    ', 'DIST_N  ', 'CTOR    ']

    ! Carte pour les fonctions
    cartcf = nomu//'.CARCOQUF'
    call exisd('CARTE', cartcf, iret)
    lcartf = iret .ne. 0

    ! Si la carte fonction n'existe pas encore, doit-on la creer ?
    if (.not. lcartf) then
        do ioc = 1, nbocc
            call getvid('GRILLE', 'SECTION_FO', iocc=ioc, scal=slf, nbret=mc_section_fo)
            call getvid('GRILLE', 'EXCENTREMENT_FO', iocc=ioc, scal=ezf, nbret=mc_excentrement_fo)
            if (mc_section_fo .ne. 0 .or. mc_excentrement_fo .ne. 0) then
                lcartf = .true.
                call alcart('V', cartcf, noma, 'CACOQU_F')
                exit
            end if
        end do
    end if

    ! Objets intermédiaires pour l'allocation de coor_allurs
    if (lcartf) then
        call jeveuo(cartcf//'.NCMP', 'E', vk8=nom_cmp_f)
        call jeveuo(cartcf//'.VALV', 'E', vk8=val_cmp_f)
        ASSERT(size(nom_cmp_f) .ge. 2)
        nom_cmp_f(1:2) = ['SECT_L  ', 'DIST_N  ']
    end if

! --------------------------------------------------------------------------------------------------
!  LECTURE DES VALEURS ET AFFECTATION DANS : CARTGR, CARTCF
! --------------------------------------------------------------------------------------------------
    do ioc = 1, nbocc

        ! Valeurs vides
        val_cmp_r(1:5) = 0.d0
        if (lcartf) val_cmp_f(1:2) = '&&ACEAGB'

        ! Affectation de la section
        call getvr8('GRILLE', 'SECTION', iocc=ioc, scal=sl, nbret=mc_section)
        if (mc_section .ne. 0) val_cmp_r(1) = sl

        call getvid('GRILLE', 'SECTION_FO', iocc=ioc, scal=slf, nbret=mc_section_fo)
        if (mc_section_fo .ne. 0) val_cmp_f(1) = slf

        ! Affectation de l'excentrement
        call getvr8('GRILLE', 'EXCENTREMENT', iocc=ioc, scal=ez, nbret=mc_excentrement)
        if (mc_excentrement .ne. 0) val_cmp_r(4) = ez

        call getvid('GRILLE', 'EXCENTREMENT_FO', iocc=ioc, scal=ezf, nbret=mc_excentrement_fo)
        if (mc_excentrement_fo .ne. 0) val_cmp_f(2) = ezf

        ! Affectation de coef_rigi_drz
        call getvr8('GRILLE', 'COEF_RIGI_DRZ', iocc=ioc, scal=val_cmp_r(5))

        ! Groupes de mailles (nombre puis lecture)
        call getvem(noma, 'GROUP_MA', 'GRILLE', 'GROUP_MA', ioc, 0, [k24bid], nbr_gma)
        nbr_gma = -nbr_gma
        ASSERT(nbr_gma .ge. 1)
        allocate (l_gma(nbr_gma))
        call getvem(noma, 'GROUP_MA', 'GRILLE', 'GROUP_MA', ioc, nbr_gma, l_gma, ibid)

        ! Paramètres d'affectation de la direction de la grille
        call getvr8('GRILLE', 'ANGL_REP_1', iocc=ioc, nbval=2, vect=ang1, nbret=mc_angl_rep_1)
        call getvr8('GRILLE', 'VECT_1', iocc=ioc, nbval=3, vect=axex, nbret=mc_vect_1)
        call getvr8('GRILLE', 'ANGL_REP_2', iocc=ioc, nbval=2, vect=ang2, nbret=mc_angl_rep_2)
        call getvr8('GRILLE', 'VECT_2', iocc=ioc, nbval=3, vect=axey, nbret=mc_vect_2)

        if (mc_vect_1 .ne. 0 .or. mc_vect_2 .ne. 0) then
            if (mc_vect_1 .ne. 0) then
                norm = sqrt(dot_product(axex, axex))
                if (norm .lt. epsi) call utmess('F', 'MODELISA_10')
                axex = axex/norm
            end if
            if (mc_vect_2 .ne. 0) then
                norm = sqrt(dot_product(axey, axey))
                if (norm .lt. epsi) call utmess('F', 'MODELISA_10')
                axey = axey/norm
            end if

            call getvtx('GRILLE', 'REPERE', iocc=ioc, scal=repere)
            if (repere .eq. 'CYLINDRIQUE') then
                call getvr8('GRILLE', 'AXE_Z', iocc=ioc, nbval=3, vect=vect_ez, nbret=mc_axe_z)
                norm = sqrt(dot_product(vect_ez, vect_ez))
                if (norm .le. epsi) call utmess('F', 'MODELISA_10')
                vect_ez = vect_ez/norm

                call getvr8('GRILLE', 'ORIGINE', iocc=ioc, nbval=3, vect=orig, nbret=mc_orig)
            end if
        end if

        ! Affectation par parametres constants (i.e. identiques pour toutes les mailles)
        if (mc_angl_rep_1 .ne. 0 .or. (mc_vect_1 .ne. 0 .and. repere .eq. 'GLOBAL')) then

            if (mc_vect_1 .ne. 0) then
                call angvx(axex, angrd(1), angrd(2))
                ang1 = angrd*r8rddg()
            end if
            val_cmp_r(2:3) = ang1(1:2)

            do igr = 1, nbr_gma
                call nocart(cartgr, 2, 5, groupma=l_gma(igr))
                if (lcartf) call nocart(cartcf, 2, 2, groupma=l_gma(igr))
            end do

            ! Affectation par parametres variables (i.e. coor_allurs différentes en chaque maille)
        else

            ! ANGL_REP_2 -> On se ramène au cas d'un vecteur directeur VECT_2
            if (mc_angl_rep_2 .ne. 0) then
                angrd = ang2*r8dgrd()
                axey(1) = cos(angrd(1))*cos(angrd(2))
                axey(2) = sin(angrd(1))*cos(angrd(2))
                axey(3) = sin(angrd(2))
                mc_angl_rep_2 = 0
                mc_vect_2 = 3
            end if

            ! Parcours des groupes de mailles
            do igr = 1, nbr_gma
                call jeveuo(jexnom(noma//'.GROUPEMA', l_gma(igr)), 'L', vi=gma)

                ! Construction des angles directeurs pour chaque maille du groupe de mailles courant
                do ima = 1, size(gma)

                    ! Type de maille
                    numa = gma(ima)
                    call jeveuo(jexnum(noma//'.CONNEX', numa), 'L', vi=l_nds)
                    tria = size(l_nds) .eq. 3 .or. size(l_nds) .eq. 6
                    quad = size(l_nds) .eq. 4 .or. size(l_nds) .eq. 8 .or. size(l_nds) .eq. 9
                    ASSERT(tria .or. quad)
                    if (tria) nno = 3
                    if (quad) nno = 4

                    ! Coordonnees des noeuds de la maille
                    do n = 1, nno
                        noe = l_nds(n)
                        coor(:, n) = coor_all(1+3*(noe-1)+1-1:1+3*(noe-1)+3-1)
                    end do

                    !  barycentre
                    bary(:) = sum(coor(:, 1:nno), dim=2)/nno

                    ! Vecteurs tangents à la maille calculés au barycentre (fonctions formes P1)
                    if (tria) then
                        vtg1 = coor(:, 2)-coor(:, 1)
                        vtg2 = coor(:, 3)-coor(:, 1)
                    else
                        vtg1 = coor(:, 2)-coor(:, 1)+coor(:, 3)-coor(:, 4)
                        vtg2 = coor(:, 3)-coor(:, 2)+coor(:, 4)-coor(:, 1)
                    end if

                    ! Vecteur normal à la maille calculé au barycentre de la maille
                    vec_n = prod_vect(vtg1, vtg2)
                    norm = sqrt(dot_product(vec_n, vec_n))
                    if (norm .le. epsi) call utmess('F', 'MODELISA_12')
                    vec_n = vec_n/norm

                    ! Passage du vecteur directeur dans le repere local si necessaire
                    if (repere .eq. 'CYLINDRIQUE') then

                        ! Repere local
                        proj = orig+dot_product(bary-orig, vect_ez)*vect_ez
                        vect_er = bary-proj
                        norm = sqrt(dot_product(vect_er, vect_er))
                        if (norm .le. epsi) call utmess('F', 'MODELISA_13')
                        vect_er = vect_er/norm
                        vect_et = prod_vect(vect_ez, vect_er)

                        ! Changement de repere
                        if (mc_vect_1 .ne. 0) vec_x = &
                            axex(1)*vect_er+axex(2)*vect_ez+axex(3)*vect_et
                        if (mc_vect_2 .ne. 0) vec_y = &
                            axey(1)*vect_er+axey(2)*vect_ez+axey(3)*vect_et
                    else
                        if (mc_vect_1 .ne. 0) vec_x = axex
                        if (mc_vect_2 .ne. 0) vec_y = axey
                    end if

                    ! Construction du vecteur directeur vec_x si nécessaire
                    if (mc_vect_2 .ne. 0) vec_x = prod_vect(vec_y, vec_n)
                    norm = sqrt(dot_product(vec_x, vec_x))
                    if (norm .le. epsi) call utmess('F', 'MODELISA_11')
                    vec_x = vec_x/norm

                    ! Calcul des angles
                    call angvx(vec_x, angrd(1), angrd(2))
                    ang1 = angrd*r8rddg()
                    val_cmp_r(2:3) = ang1

                    ! Affectation dans la carte
                    call nocart(cartgr, 3, 5, mode='NUM', nma=1, limanu=[numa])
                    if (lcartf) call nocart(cartcf, 3, 2, mode='NUM', nma=1, limanu=[numa])
                end do
            end do
        end if
        deallocate (l_gma)
    end do

    ! Si pas membrane
    if (.not. locamb) then
        call jedetr(cartgr//'.NCMP')
        call jedetr(cartgr//'.VALV')
        if (lcartf) then
            call jedetr(cartcf//'.NCMP')
            call jedetr(cartcf//'.VALV')
        end if
    end if

    call jedema()
end subroutine
