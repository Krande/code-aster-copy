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
subroutine utmasu(mail, kdim, nbCell, listCellNume, nomob1, &
                  coor, nbmavo, mailvo, coince)
    implicit none
! person_in_charge: jacques.pellet at edf.fr
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/assert.h"
#include "asterfort/indiis.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/oriem0.h"
#include "asterfort/utmavo.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/int_to_char8.h"
!
    integer(kind=8) :: nbCell, nbmavo, mailvo(*)
    real(kind=8) :: coor(*)
    integer(kind=8), pointer :: listCellNume(:)
    character(len=2) :: kdim
    character(len=8) :: mail
    character(len=*) :: nomob1
    aster_logical :: coince
! but : recuperer la couche des mailles qui "bordent"
!       une liste de mailles de "peau"
!
!     peau 2d (tria,quad) => mailles 3d
!     peau 1d (seg)       => mailles 2d
!
!   arguments :
!   -----------
!     mail (in)  : nom du maillage
!     kdim (in)  : / '3D' recherche les mailles 3d voisines
!                  / '2D' recherche les mailles 2d voisines
!     nbCell (in)  : nombre de mailles de lima
!     lima  (in)  : liste des numeros des mailles de peau
!
!     nomob1 (in/jxout) : nom de l' ojb a creer (vecteur d'entiers)
!       Ce vecteur est de longueur nbCell.
!       Pour chaque maille de peau, il contient une maille qui
!       "borde" cette maille de peau.
!       0 : si la maille de peau n'a pas de maille qui la borde.
!
!     nbmavo (in) : / nb de mailles de mailvo
!                   / 0
!     mailvo (in) : / liste de maille "candidates"
!                   / 0 (si nbmavo=0)
!     mailvo est un sous-ensemble des mailles du maillage
!     qui permet d'"orienter" le choix des mailles de nomob1
!     mailvo est utilise en particulier quand on veut reorienter
!     des facettes qui sont inserrees entre des mailles volumiques:
!        orie_peau_3d ("group_ma_volu")
!        orie_peau_2d ("group_ma_surf")
!
!     coince (in) : /.true.  /.false.
!       Si .true. on accepte qu'une maille de peau soit "coincee"
!       entre 2 mailles de part et d'autre de la peau.
!       Sinon, on emet une erreur fatale.
!       si .true._1, on choisit la maille telle que sa normale sortante
!       soit la meme que celle de la maille de peau.
!
!-----------------------------------------------------------------------
!
    integer(kind=8) :: p1, p2, p3, p4, jm3d, nbmat, im1, im2
    integer(kind=8) :: iCell, cellNume, nnoe, ino, nbm, i, k, indi, nnoem, nnoe1
    integer(kind=8) :: ifm, niv, ipos, cellTypeNume
    integer(kind=8) :: lisnoe(27), indmai
    aster_logical :: first
    character(len=8) :: k8b, nomail, cellTypeName
    character(len=16) :: oper, k16b
    character(len=24) :: nomavo, valk(4)
    integer(kind=8), pointer :: typmail(:) => null()
!     ------------------------------------------------------------------
    call jemarq()
    call infniv(ifm, niv)
!
    ASSERT(kdim .eq. '3D' .or. kdim .eq. '2D')
    first = .false.
    call getres(k8b, k16b, oper)
!
    call jeveuo(jexatr(mail//'.CONNEX', 'LONCUM'), 'L', p2)
    call jeveuo(mail//'.CONNEX', 'L', p1)
    call jeveuo(mail//'.TYPMAIL', 'L', vi=typmail)
!
! --- CREATION DE NOMOB1:
!     -------------------
    call wkvect(nomob1, 'V V I', nbCell, jm3d)
!
!   -- recuperation des mailles voisines de lima :
!   ----------------------------------------------
    nomavo = '&&UTMASU.MAILLE_VOISINE '
    call utmavo(mail, kdim, listCellNume, nbCell, 'V', &
                nomavo, nbmavo, mailvo)
    call jeveuo(jexatr(nomavo, 'LONCUM'), 'L', p4)
    call jeveuo(nomavo, 'L', p3)
!
!
!   -- on remplit nomob1 :
!   ----------------------
    do iCell = 1, nbCell
        cellNume = listCellNume(iCell)
        cellTypeNume = typmail(cellNume)
        nnoe = zi(p2+cellNume)-zi(p2-1+cellNume)
        ASSERT(nnoe .le. 27)
        do ino = 1, nnoe
            lisnoe(ino) = zi(p1-1+zi(p2+cellNume-1)+ino-1)
        end do
        nbmat = zi(p4+iCell+1-1)-zi(p4+iCell-1)
        nbm = 0
        do i = 1, nbmat
            im2 = zi(p3+zi(p4+iCell-1)-1+i-1)
            if (im2 .eq. 0) goto 10
            if (zi(p1+zi(p2+im2-1)-1) .eq. 0) goto 10
            nnoem = zi(p2+im2)-zi(p2-1+im2)
!
            do k = 1, nnoe
                indi = indiis(zi(p1+zi(p2+im2-1)-1), lisnoe(k), 1, nnoem)
                if (indi .eq. 0) goto 10
            end do
            nbm = nbm+1
            if (nbm .eq. 1) then
                zi(jm3d+iCell-1) = im2
            else
!               -- cas ou la maille de peau est bordee par plus
!                  d'une maille. Il faut verifier la coherence.
                im1 = zi(jm3d+iCell-1)
                nnoe1 = zi(p2+im1)-zi(p2-1+im1)
                call jenuno(jexnum('&CATA.TM.NOMTM', cellTypeNume), cellTypeName)
                call oriem0(kdim, cellTypeName, coor, zi(p1+zi(p2+im1-1)-1), nnoe1, &
                            zi(p1+zi(p2+im2-1)-1), nnoem, lisnoe, nnoe, ipos, &
                            indmai)
                if (ipos .eq. 0) then
!                   -- si ipos=0, les 2 mailles sont du meme cote, on conserve la 1ere.
                else
                    if (indmai .lt. 0) then
                        ASSERT(indmai .eq. -1 .or. indmai .eq. -2 .or. indmai .eq. -12)
                        if (indmai .eq. -1) then
                            valk(1) = int_to_char8(im1)
                            call utmess('F', 'CALCULEL2_32', sk=valk(1))
!                           -- C'est la maille 1 qui est degeneree => on prend la 2 :
                            zi(jm3d+iCell-1) = im2
                        else
                            valk(1) = int_to_char8(im2)
                            call utmess('F', 'CALCULEL2_32', sk=valk(1))
                        end if
!
                    else
!
!                       -- sinon, im2 et im1 sont de part et d'autre de cellNume
                        if (.not. coince) then
                            valk(1) = int_to_char8(cellNume)
                            valk(2) = int_to_char8(im1)
                            valk(3) = int_to_char8(im2)
                            call utmess('F', 'PREPOST4_97', nk=3, valk=valk)
                        else
                            ASSERT(indmai .eq. 1 .or. indmai .eq. 2)
                            if (indmai .eq. 2) then
!                               -- C'est im2 (qui est du cote "-" de iCell) qu'il faut retenir :
                                zi(jm3d+iCell-1) = im2
                            end if
                        end if
                    end if
                end if
            end if
10          continue
        end do
!
!
        if (nbm .eq. 0 .and. niv .gt. 1) then
            nomail = int_to_char8(cellNume)
            if (first) then
                valk(1) = nomail
                call utmess('A+', 'PREPOST6_29', sk=valk(1))
            else
                valk(1) = nomail
                call utmess('A+', 'PREPOST6_30', sk=valk(1))
            end if
            first = .true.
        end if
!
    end do
!
    call jedetr(nomavo)
    call jedema()
!
end subroutine
