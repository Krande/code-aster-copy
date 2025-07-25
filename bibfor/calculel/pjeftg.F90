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

subroutine pjeftg(igeom, geomi, nomai, motfac, iocc, ncas)
!
!
! --------------------------------------------------------------------------------------------------
!
!   Transformer la géometrie des noeuds du MAILLAGE_2 avant la la projection.
!       Cela permet par exemple de projeter :
!           - un mailllage sur un autre maillage homothetique
!           - un maillage 2d sur un maillage 3d "ecrasé" dans un plan (2d axis -> 3d axis)
!
!   in :
!       igeom       : numéro de la géométrie à transformer 1 ou 2
!       nomai       : nom du maillage initial
!
!   out :
!       geomi (k24) : nom de l'objet contenant la géométrie transformée du MAILLAGE_[1|2]
!                     si pas de géométrie transformée geomi=' '
!
! --------------------------------------------------------------------------------------------------
!
    implicit none
    integer(kind=8) :: igeom
    character(len=8) :: nomai
    character(len=24) :: geomi
    character(len=*) :: motfac
    integer(kind=8) :: iocc
    character(len=8) :: ncas
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/assert.h"
#include "asterfort/copisd.h"
#include "asterfort/fointe.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/infniv.h"
#include "asterfort/irmail.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedupo.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lxlgut.h"
#include "asterfort/ulaffe.h"
#include "asterfort/ulnume.h"
#include "asterfort/ulopen.h"
#include "asterfort/utmess.h"
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: n1, nbnoi, ifonc, ibid, niveau, nkcar
    integer(kind=8) :: nfonc, jgeomi, inoi, ier, unite
    real(kind=8) :: vx(3)
    character(len=1) :: kgeom
    character(len=8) :: lfonc(3), lparx(3), maili, k8bid
    character(len=16) :: formar, method
    character(len=16) :: concept, nomcmd
    character(len=19) :: resout
    character(len=80) :: fichier
    aster_logical :: limpr
!
    character(len=14) :: messk(3)
!
    character(len=14) :: trans12(2), geome12(2)
    data trans12/'TRANSF_GEOM_1', 'TRANSF_GEOM_2'/
    data geome12/'&&PJEFTG.GEOM1', '&&PJEFTG.GEOM2'/
!
! --------------------------------------------------------------------------------------------------
    call jemarq()
!
    ASSERT((igeom .eq. 1) .or. (igeom .eq. 2))
!
!   Si motfac<>' ' On vérifie que TRANSF_GEOM_* n'a pas été utilisé en mot clef simple "au dessus"
    if (motfac .ne. ' ') then
        call getvid(' ', trans12(igeom), nbret=nfonc)
        if (nfonc .ne. 0) then
            messk(1) = motfac
            messk(2) = trans12(igeom)
            call utmess('F', 'PROJECTION4_3', nk=2, valk=messk)
        end if
    end if
!
!   Prise en compte du mot-cle TRANSF_GEOM_[1|2]
    call getvid(motfac, trans12(igeom), iocc=iocc, nbval=3, vect=lfonc, nbret=nfonc)
    ASSERT(nfonc .ge. 0)
!   Transformation de la géométrie ou pas
    geomi = ' '
    maili = ' '
    if (nfonc .gt. 0) then
        ASSERT(nfonc .eq. 2 .or. nfonc .eq. 3)
        if (nfonc .eq. 2) then
            if (ncas .eq. '2D') then
                lfonc(3) = '&FOZERO'
            else
                messk(1) = motfac
                messk(2) = trans12(igeom)
                messk(3) = ncas
                call utmess('F', 'PROJECTION4_4', nk=3, valk=messk)
            end if
        end if
        geomi = geome12(igeom)
        call jedetr(geomi)
!       Copiage du maillage initial
        maili = '&&PJEFTG'
        call jedetr(maili)
        call copisd('MAILLAGE', 'V', nomai, maili)
!       Nombre de noeuds, adresse des coordonnées
        call jelira(maili//'.COORDO    .VALE', 'LONMAX', n1)
        call jeveuo(maili//'.COORDO    .VALE', 'E', jgeomi)
        nbnoi = n1/3
        ASSERT(n1 .eq. nbnoi*3)
        lparx(1) = 'X'
        lparx(2) = 'Y'
        lparx(3) = 'Z'
!       Modification des coordonnées de 'maili'
        do inoi = 1, nbnoi
            do ifonc = 1, 3
                call fointe('F', lfonc(ifonc), 3, lparx, zr(jgeomi+3*(inoi-1)), vx(ifonc), ier)
                ASSERT(ier .eq. 0)
            end do
            zr(jgeomi-1+3*(inoi-1)+1) = vx(1)
            zr(jgeomi-1+3*(inoi-1)+2) = vx(2)
            zr(jgeomi-1+3*(inoi-1)+3) = vx(3)
        end do
!       Copiage de la géométrie modifiée dans geomi
        call jedupo(maili//'.COORDO    .VALE', 'V', geomi, ASTER_FALSE)
    end if

!   Si INFO=2, on imprime au format MED les maillages utilisés lors de la projection
!   si ceux-ci ne sont pas connus de l'utilisateur, c'est à dire si :
!     METHODE='SOUS_POINT' (pour le maillage "2")
!     et/ou si l'utilisateur a demandé TRANSF_GEOM_1/_2
    call infniv(ibid, niveau)
    if (niveau .le. 1) goto 999

    call getvtx(' ', 'METHODE', scal=method, nbret=ibid)
    if (ibid .eq. 0) method = ' '

    limpr = .false.
    if (method .eq. 'SOUS_POINT') then
        if ((nfonc .gt. 0) .or. (igeom .eq. 2)) limpr = .true.
    else
        if (nfonc .gt. 0) limpr = .true.
    end if

    if (limpr) then
        ibid = 0
        k8bid = '        '
        formar = '        '
!       Réservation d'une unité et codage du nom en dur
!           Récupération du nom de la commande
!               le nom utilisateur du résultat : resout
!               le nom du concept résultat     : concept (en cas de reuse)
!               le nom de la commande          : nomcmd
        call getres(resout, concept, nomcmd)
!       numéro de la géométrie
        write (kgeom, '(I1)') igeom
        nkcar = lxlgut(resout)
        fichier = 'REPE_OUT/'//resout(1:nkcar)//'_'//kgeom//'.med'
!       unite libre
        unite = ulnume()
        if (unite .le. 0) call utmess('F', 'UTILITAI5_10')
!       On ouvre et écriture
        call ulaffe(unite, fichier, ' ', 'N', 'O')
        if (geomi .ne. ' ') then
            call irmail('MED', unite, ibid, maili, ASTER_FALSE, k8bid, 1, formar)
        else
            call irmail('MED', unite, ibid, nomai, ASTER_FALSE, k8bid, 1, formar)
        end if
!       On ferme
        call ulopen(-unite, k8bid, k8bid, k8bid, k8bid)
    end if
!
999 continue
    if (maili .ne. ' ') call jedetr(maili)
    call jedema()
end subroutine
