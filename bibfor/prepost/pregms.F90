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
subroutine pregms(igmsh, imod)
! aslint: disable=
    implicit none
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/codent.h"
#include "asterfort/gmeelt.h"
#include "asterfort/gmeneu.h"
#include "asterfort/gmlelt.h"
#include "asterfort/gmlneu.h"
#include "asterfort/inigms.h"
#include "asterfort/iunifi.h"
#include "asterfort/jedetr.h"
#include "asterfort/jjmmaa.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"

    integer(kind=8) :: igmsh, imod
!.======================================================================
!
!      PREGMS --   INTERFACE GMSH --> ASTER
!                  LECTURE DU FICHIER  .GMSH
!                  ECRITURE DU FICHIER .MAIL
!
!   ARGUMENT        E/S  TYPE         ROLE
!    IGMSH          IN    I         UNITE LOGIQUE DU FICHIER GMSH
!    IMOD           IN    I         UNITE LOGIQUE DU FICHIER MAIL
!
! ......................................................................
!
!
!
!
    character(len=4) :: ct(3)
    character(len=8) :: rquoi
    character(len=12) :: aut, debfic, finnod, debelm
    character(len=14) :: aut1, current_line
    integer(kind=8) :: i, imes, nbmail, nbnode, versio, maxnod, nbtyma
    integer(kind=8) :: vali(1), nbgrou, gr_tag, ibid
    aster_logical :: exigr
!
    parameter(maxnod=32, nbtyma=19)
    integer(kind=8) :: nbnoma(nbtyma), nuconn(nbtyma, maxnod), i_gr, ipos
    integer(kind=8) :: gr_dim, igr, ima
    character(len=8) :: nomail(nbtyma), gr_name, gr_name_raw, group_name
!
    integer(kind=8), pointer :: list_nums_gr(:) => null()
    integer(kind=8), pointer :: list_dims_gr(:) => null()
    integer(kind=8), pointer :: list_indice_gr(:) => null()
    character(len=8), pointer :: list_noms_gr(:) => null()
!
! ----------------------------------------------------------------------
!8
! ---- INITIALISATIONS
!      ---------------
    rquoi = '????????'
!
    do i = 1, nbtyma
        nomail(i) = rquoi
    end do
!
! --- RECUPERATION DES NUMEROS D'UNITE LOGIQUE DES FICHIERS :
!     -----------------------------------------------------
    imes = iunifi('MESSAGE')
!
! --- AFFECTATION DE NOMAIL AVEC LE NOM DU TYPE DES ELEMENTS :
!     ------------------------------------------------------
    call inigms(nomail, nbnoma, nuconn)
!
! --- LECTURE EN DEBUT DU FICHIER .GMSH POUR DETERMINER LE FORMAT :
!     -----------------------------------------
    read (igmsh, *) debfic
!
    if (debfic(1:4) .eq. '$NOD') then
        versio = 1
        !
        ! Read the one group to which each element belongs
        read (igmsh, *) current_line
        do while (current_line(1:4) /= '$ELM')
            read (igmsh, *) current_line
        end do
        read (igmsh, '(I10)') nbmail
        call wkvect('&&PREGMS.INDICE.GROUP_MA', 'V V I', nbmail, &
                    vi=list_indice_gr)
        do ima = 1, nbmail
            read (igmsh, *) ibid, ibid, list_indice_gr(ima)
        end do
        !
        ! Determine number of groups
        nbgrou = 1
        do ima = 2, nbmail
            exigr = .false._1
            do i = 1, ima-1
                if (list_indice_gr(i) == list_indice_gr(ima)) then
                    exigr = .true._1
                    exit
                end if
            end do
            if (.not. exigr) nbgrou = nbgrou+1
        end do
        !
        ! Fill in ids and (automatic) names of groups
        call wkvect('&&PREGMS.NUMERO.GROUP_MA', 'V V I', nbgrou, &
                    vi=list_nums_gr)
        call wkvect('&&PREGMS.NOMS.GROUP_MA', 'V V K8', nbgrou, &
                    vk8=list_noms_gr)
        list_nums_gr(1) = list_indice_gr(1)
        group_name(1:2) = 'GM'
        call codent(list_indice_gr(1), 'G', group_name(3:8))
        list_noms_gr(1) = group_name
        igr = 1
        do ima = 2, nbmail
            exigr = .false._1
            do i = 1, ima-1
                if (list_indice_gr(i) == list_indice_gr(ima)) then
                    exigr = .true._1
                    exit
                end if
            end do
            if (.not. exigr) then
                igr = igr+1
                list_nums_gr(igr) = list_indice_gr(ima)
                group_name(1:2) = 'GM'
                call codent(list_indice_gr(ima), 'G', group_name(3:8))
                list_noms_gr(igr) = group_name
            end if
        end do
        !
        ! Comme une ancienne cassette
        rewind igmsh
        read (igmsh, *) current_line

    else if (debfic(1:11) .eq. '$MeshFormat') then
        ! Dont forget to initialize nbgrou in case there is none
        nbgrou = 0
        versio = 2
        read (igmsh, *) current_line
        do while (current_line(1:6) /= '$Nodes' .and. &
                  current_line(1:14) /= '$PhysicalNames')
            read (igmsh, *) current_line
        end do
        if (current_line(1:14) == '$PhysicalNames') then
            read (igmsh, *) nbgrou
            write (imes, *)
            write (imes, *) 'LECTURE DES NOMS DES GROUPES'
            !
            call wkvect('&&PREGMS.NUMERO.GROUP_MA', 'V V I', nbgrou, &
                        vi=list_nums_gr)
            call wkvect('&&PREGMS.DIME.GROUP_MA', 'V V I', nbgrou, &
                        vi=list_dims_gr)
            call wkvect('&&PREGMS.NOMS.GROUP_MA', 'V V K8', nbgrou, &
                        vk8=list_noms_gr)
            do i_gr = 1, nbgrou
                read (igmsh, *) gr_dim, gr_tag, gr_name_raw
                ipos = index(gr_name_raw, '"', .true.)
                if (ipos == 0) gr_name = gr_name_raw(1:8)
                if (ipos /= 0) gr_name = gr_name_raw(1:ipos-1)
                list_dims_gr(i_gr) = gr_dim
                list_nums_gr(i_gr) = gr_tag
                list_noms_gr(i_gr) = gr_name
            end do
            do while (current_line(1:6) /= '$Nodes')
                read (igmsh, *) current_line
            end do
        end if
    else
        call utmess('F', 'PREPOST6_38')
    end if
!
    call utmess('I', 'PREPOST6_39')
    vali(1) = versio
    call utmess('I', 'PREPOST6_40', si=vali(1))
!
! --- ECRITURE DU TITRE DANS LE FICHIER .MAIL :
!     ---------------------------------------
    write (imod, '(A)') 'TITRE'
!
! --- ECRITURE DE LA DATE DU JOUR :
!     ---------------------------
    call jjmmaa(ct, aut)
    aut1 = 'INTERFACE_GMSH'
    write (imod, '(9X,2A,17X,A,A2,A,A2,A,A4)') 'AUTEUR=', aut1, 'DATE=',&
     &  ct(1) (1:2), '/', ct(2) (1:2), '/', ct(3)
    write (imod, '(A)') 'FINSF'
    write (imod, '(A)') '%'
!
! --- LECTURE DES NOEUDS ET DE LEURS COORDONNEES DANS LE FICHIER .GMSH:
!     ----------------------------------------------------------------
    write (imes, *)
    write (imes, *) 'LECTURE DES NOEUDS ET DE LEURS COORDONNEES'
    call gmlneu(igmsh, nbnode)
!
! --- FIN DE LA LECTURE DES NOEUDS :
!     ----------------------------
    read (igmsh, *) finnod
!
    if ((finnod(1:7) .ne. '$ENDNOD') .and. (finnod(1:9) .ne. '$EndNodes')) then
        call utmess('F', 'PREPOST6_41')
    end if
!
! --- DEBUT DE LA LECTURE DES ELEMENTS DANS LE FICHIER .GMSH :
!     ------------------------------------------------------
    write (imes, *)
    write (imes, *) 'LECTURE DES MAILLES'
    read (igmsh, *) debelm
!
    if ((debelm(1:4) .ne. '$ELM') .and. (debelm(1:9) .ne. '$Elements')) then
        call utmess('F', 'PREPOST6_42')
    end if
!
!
! --- LECTURE DES MAILLES ET DES GROUP_MA :
!     -----------------------------------
    call gmlelt(igmsh, maxnod, nbtyma, nbmail, nbnoma, &
                nuconn, versio, nbgrou)
!
! --- ECRITURE DES NOEUDS ET DE LEURS COORDONNEES DANS LE FICHIER .MAIL:
!     -----------------------------------------------------------------
    call gmeneu(imod, nbnode)
!
! --- ECRITURE DES MAILLES ET DES GROUP_MA DANS LE FICHIER .MAIL :
!     ----------------------------------------------------------
    call gmeelt(imod, nbtyma, nomail, nbnoma, nuconn, &
                nbmail, nbgrou)
!
! --- MENAGE :
!     ------
    call jedetr('&&PREGMS.INFO.NOEUDS')
    call jedetr('&&PREGMS.COOR.NOEUDS')
    call jedetr('&&PREGMS.NUMERO.MAILLES')
    call jedetr('&&PREGMS.TYPE.MAILLES')
    call jedetr('&&PREGMS.GROUPE.MAILLES')
    call jedetr('&&PREGMS.NBNO.MAILLES')
    call jedetr('&&PREGMS.CONNEC.MAILLES')
    call jedetr('&&PREGMS.NBMA.GROUP_MA')
    call jedetr('&&PREGMS.NBTYP.MAILLES')
    call jedetr('&&PREGMS.LISTE.GROUP_MA')
    call jedetr('&&PREGMS.INDICE.GROUP_MA')
    call jedetr('&&PREGMS.GRMA.MAILLES')
    if (versio == 2) then
        call jedetr('&&PREGMS.NOMS.GROUP_MA')
        call jedetr('&&PREGMS.DIME.GROUP_MA')
        call jedetr('&&PREGMS.NUMERO.GROUP_MA')
    end if
!
!.============================ FIN DE LA ROUTINE ======================
end subroutine
