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
subroutine fonmai2(resu, nomail, typfon, iocc, nbnoff, &
                   typm)
    implicit none
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/cgnoor.h"
#include "asterfort/getvtx.h"
#include "asterfort/i2extf.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/ornofd.h"
#include "asterfort/utmess.h"
#include "asterfort/utnono.h"
#include "asterfort/wkvect.h"
#include "asterfort/char8_to_int.h"
#include "asterfort/int_to_char8.h"
!
    integer(kind=8) :: iocc, nbnoff
    character(len=8) :: resu, nomail, typfon
! FONCTION REALISEE:
!
!     VERIFICATION DES ENTITES LORSQUE LE FRONT EST DECRIT PAR
!     DES GROUPES DE MAILLES RENSEIGNEES DANS DEFI_FOND_FISS
!     CONSTRUCTION DU FOND DE FISSURE A PARTIR CES DONNEES
!
!     ENTREES:
!        RESU   : NOM DU CONCEPT RESULTAT DE L'OPERATEUR
!        NOMAIL : NOM DU MAILLAGE
!        TYPFON : TYPE DE FOND
!                 IL PEUT VALOIR OUVERT/FERME
!        IOCC   : OCCURENCE COURANTE DE MOTFAC
!     SORTIES:
!        NBNOFF : NOMBRE DE NOEUDS EN FOND DE FISSURE
!        TYPM   : FOND LINEAIRE OU QUADRATIQUE : SEG2 OU SEG3
!
!-----------------------------------------------------------------------
!
    real(kind=8) :: vecori(3)
!
    integer(kind=8) :: jcour2, iatyma, idnono, idlino
    integer(kind=8) :: i, nbma, n1, im, nig
    integer(kind=8) :: nid, numno, iret, trouv, numma
    character(len=8) :: k8b, nomma, typm, ndorig, ndextr
    character(len=8) :: noeud
    character(len=16) :: k16bid, nomcmd, motfac
    character(len=16) :: motcle(2), typmcl(2)
    character(len=24) :: conec, typp, noeord
    character(len=24) :: mesnoe, mafour, nogrp
    integer(kind=8), pointer :: maillestriees(:) => null()
! DEB-------------------------------------------------------------------
    call jemarq()
!
    call getres(k8b, k16bid, nomcmd)
!     ------------------------------------------------------------------
!     INITIALISATION DE VARIABLES
!     ------------------------------------------------------------------
    motfac = 'FOND_FISS'
    typp = nomail//'.TYPMAIL        '
    conec = nomail//'.CONNEX         '
    call jeveuo(typp, 'L', iatyma)
    call jenuno(jexnum('&CATA.TM.NOMTM', zi(iatyma)), typm)
!
!     ------------------------------------------------------------------
!     --- RECHERCHE DES NOEUDS SOMMET DES MAILLES RENSEIGNEES
!     ------------------------------------------------------------------
    motcle(1) = 'GROUP_MA'
    typmcl(1) = 'GROUP_MA'
    mafour = '&&FONMAI.MALIGNE'
    call cgnoor(mafour, nomail, motfac, iocc, 1, &
                motcle, typmcl, typfon, nbma, ndorig, &
                ndextr, typm, vecori)
    call jeveuo(mafour, 'L', jcour2)
!
!
!     ------------------------------------------------------------------
!     --- SI FERME : RECUPERATION DE MAILLE_ORIG POUR AVOIR
!     --- LE SENS DE PARCOURS DE LA COURBE FERMEE
!     ------------------------------------------------------------------
!
    if (typfon .eq. 'FERME') then
!
        numma = 0
        call getvtx(motfac, 'GROUP_MA_ORIG', iocc=1, nbval=0, nbret=n1)
        if (n1 .ne. 0) then
            call getvtx(motfac, 'GROUP_MA_ORIG', iocc=1, scal=nogrp, nbret=n1)
            call utnono(' ', nomail, 'MAILLE', nogrp, nomma, &
                        iret)
            if (iret .eq. 10) then
                call utmess('F', 'RUPTURE0_41', sk=nogrp)
            else if (iret .eq. 1) then
                call utmess('F', 'RUPTURE0_45', sk=ndorig)
            end if
            numma = char8_to_int(nomma)
        end if
!
        if (numma .eq. 0) then
            call utmess('F', 'RUPTURE0_42')
        else
            numno = char8_to_int(ndorig)
            call i2extf(numma, 1, conec(1:15), typp(1:16), nig, &
                        nid)
            if ((numno .ne. nig) .and. (numno .ne. nid)) then
                call utmess('F', 'RUPTURE0_43')
            end if
            trouv = 0
            do im = 1, nbma
                if (numma .eq. zi(jcour2-1+im)) trouv = im
            end do
            if (trouv .eq. 0) then
                call utmess('F', 'RUPTURE0_44', sk=nomma)
            else
!
!     ON REMONTE LA MAILLE_ORIG EN TETE DE LISTE
!
                AS_ALLOCATE(vi=maillestriees, size=3*nbma)
                do im = trouv, nbma
                    maillestriees(im+1-trouv) = zi(jcour2-1+im)
                end do
                do im = 1, trouv-1
                    maillestriees(im+1+nbma-trouv) = zi(jcour2-1+im)
                end do
                do im = 1, nbma
                    zi(jcour2-1+im) = maillestriees(im)
                end do
                AS_DEALLOCATE(vi=maillestriees)
            end if
        end if
    end if
!
!     ------------------------------------------------------------------
!     --- ORDONNANCEMENT DES NOEUDS EN FOND DE FISSURE
!     ------------------------------------------------------------------
    mesnoe = '&&FONMAI.NOEUD'
    call ornofd(mafour, nomail, nbma, mesnoe, ndorig, &
                ndextr, 'V', vecori)
    noeord = resu//'.FOND.NOEU'
    call jelira(mesnoe, 'LONMAX', nbnoff)
    call jeveuo(mesnoe, 'L', idnono)
!
    call wkvect(noeord, 'G V K8', nbnoff, idlino)
    do i = 1, nbnoff
        noeud = int_to_char8(zi(idnono-1+i))
        zk8(idlino-1+i) = noeud
    end do
!
!     ------------------------------------------------------------------
    call jedetr(mesnoe)
    call jedetr(mafour)
!
    call jedema()
end subroutine
