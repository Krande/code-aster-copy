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

subroutine gver2d(nocc, noeud, rinf, rsup)
    implicit none
!
!     ------------------------------------------------------------------
!
! FONCTION REALISEE:
!
!     MOT CLE FACTEUR THETA:
!
!     POUR LE NOEUD DU FOND DE FISSURE ON RECUPERE
!     LE TRIPLET ( MODULE(THETA), R_INF, R_SUP )
!
!     PUIS ON VERIFIE:
!                     QUE LE NOM DU GROUPE OU D'ELEMENTS (NOEUD)
!                     APPARTIENNENT BIEN AU MAILLAGE
!
!     ------------------------------------------------------------------
! ENTREE:
!
!     NOMA   : NOM DU MAILLAGE
!     NOCC   : NOMBRE D'OCCURENCES
!     NOMNO  : NOM DE L'OBJET CONTENANT LES NOMS DES NOEUDS
!
! SORTIE:
!
!     NOEUD      : NOEUD DU FOND DE FISSURE
!     R_INF       : RAYON INFERIEUR DE LA COURONNE
!     R_SUP       : RAYON SUPERIEUR DE LA COURONNE
!
!     ------------------------------------------------------------------
!
#include "jeveux.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
!
    character(len=8) :: config, noeud, fond, kfon
    character(len=24) :: chfond, taillr
!
    integer(kind=8) :: iocc, nocc, n1
    integer(kind=8) :: nbm, n2, lnoff, numfon, ibid
    integer(kind=8) :: iatmno
!
    real(kind=8) :: rbid, rinf, rsup, valr(2)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    call jemarq()
!
    do iocc = 1, nocc
!
        call getvr8('THETA', 'R_INF', iocc=iocc, scal=rinf, nbret=nbm)
        call getvr8('THETA', 'R_SUP', iocc=iocc, scal=rsup, nbret=nbm)
!
        if (nbm .ne. 0 .and. rsup .le. rinf) then
            call utmess('F', 'RUPTURE1_6')
        end if
!
        call getvr8('THETA', 'R_INF_FO', iocc=iocc, scal=rbid, nbret=ibid)
        if (ibid .ne. 0) then
            call utmess('F', 'RUPTURE1_18')
        end if
!
        call getvid('THETA', 'FOND_FISS', iocc=1, scal=fond, nbret=n1)
        if (n1 .ne. 0) then
!           CAS CLASSIQUE
            chfond = fond//'.FOND.NOEU'
            call jelira(chfond, 'LONMAX', lnoff)
            if (lnoff .ne. 1) then
                call utmess('F', 'RUPTURE1_10')
            else
                call jeveuo(chfond, 'L', n1)
                noeud = zk8(n1)
            end if
            numfon = 1
            if (nbm .eq. 0) then
                call dismoi('CONFIG_INIT', fond, 'FOND_FISS', repk=config)
                if (config .eq. 'DECOLLEE') then
                    call utmess('F', 'RUPTURE1_7')
                end if
            end if
        else
!           CAS X-FEM
            call getvid('THETA', 'FISSURE', iocc=1, scal=fond, nbret=n2)
            if (n2 .eq. 0) then
                call utmess('F', 'RUPTURE1_11')
            end if
!           RECUPERATION DU NUMERO DU FOND DE FISSURE DEMANDE
            call getvis('THETA', 'NUME_FOND', iocc=1, scal=numfon, nbret=ibid)
!           ON ECRIT 'NUM'+_i OU i=NUMFON
!           A LA PLACE DU NOM DU NOEUD EN FOND DE FISSURE
            call codent(numfon, 'G', kfon)
            noeud(1:8) = 'NUM_'//kfon
        end if
!
!         RECUPERATION DE RINF ET DE RSUP DANS LA SD
        if (nbm .eq. 0) then
            taillr = fond//'.FOND.TAILLE_R'
            call jeveuo(taillr, 'L', iatmno)
            rinf = 2*zr(iatmno-1+numfon)
            rsup = 4*zr(iatmno-1+numfon)
            valr(1) = rinf
            valr(2) = rsup
            call utmess('I', 'RUPTURE1_5', nr=2, valr=valr)
        end if
    end do
!
    call jedema()
end subroutine
